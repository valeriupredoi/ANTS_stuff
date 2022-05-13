## Issue in ANTS

### Problem Overview

Full trace:

```
TopologyException: Input geom 0 is invalid: Self-intersection at or near point 187.71999056058198 -28.984452570559313 at 187.71999056058198 -28.984452570559313
Traceback (most recent call last):
  File "ancil_lct.py", line 274, in <module>
    args.land_threshold,
  File "ancil_lct.py", line 206, in main
    land_fraction_threshold,
  File "ancil_lct.py", line 139, in load_data
    target_cube = target_cube.extract(ants.ExtractConstraint(source, fix_period=True))
  File "/opt/conda/envs/ants/lib/python3.7/site-packages/iris/cube.py", line 2306, in extract
    return constraint.extract(self)
  File "/opt/conda/envs/ants/lib/python3.7/site-packages/ants/_constraints.py", line 201, in extract
    source, self._target, fix_period=self._fix_period, pad_width=self._pad_width
  File "/opt/conda/envs/ants/lib/python3.7/site-packages/ants/_constraints.py", line 113, in _extract_overlap
    box = _bounding_box(target_x, target_y, src_crs)
  File "/opt/conda/envs/ants/lib/python3.7/site-packages/ants/_constraints.py", line 88, in _bounding_box
    boxes = unary_union(boxes)
  File "/opt/conda/envs/ants/lib/python3.7/site-packages/shapely/ops.py", line 149, in unary_union
    return geom_factory(lgeos.methods['unary_union'](collection))
  File "/opt/conda/envs/ants/lib/python3.7/site-packages/shapely/geometry/base.py", line 76, in geom_factory
    raise ValueError("No Shapely geometry can be created from null value")
ValueError: No Shapely geometry can be created from null value
[FAIL] python_ants -s ancil_lct.py ${source} \
[FAIL] --target-grid ${target_grid} --transform-path ${transformpath} \
[FAIL] -o ${output_vegfrac} --landseamask-output ${output_lsm}       \
[FAIL] --ants-config ${ANTS_CONFIG} <<'__STDIN__'
[FAIL] 
[FAIL] '__STDIN__' # return-code=1
2022-05-12T09:08:36Z CRITICAL - failed/EXIT
```

Source of problem:

In `ants/_constraints.py: _extract_overlap(source, target, fix_period=False, pad_width=0)._bounding_box(target_x, target_y, src_crs): l.25`:

```python
            boxes = sorted(boxes, key=lambda box: box.bounds[0])
            if getattr(target_x.units, "modulus", None):
                if utils.ndarray.allclose(
                    boxes[0].bounds[0] % target_x.units.modulus,
                    boxes[1].bounds[2] % target_x.units.modulus,
                ):
                    boxes[0] = shapely.affinity.translate(
                        boxes[0], xoff=target_x.units.modulus
                    )
                    boxes = unary_union(boxes)
```

None of the `box` in `boxes` is valid ie `box.is_valid` returns `False`.

However, there is no check for that, and the code fails with a nondescriptive error from `shapely`.

### Suggested fix

Catch the instance of any of the boxes being an invalid Shapely geometry, and fix it if you want to,
or raise a descriptive exception, e.g. a fix on-the-fly applies the fix in this [StackOverFlow Issue](https://stackoverflow.com/questions/31391209/valueerror-no-shapely-geometry-can-be-created-from-null-value)

```python
        if len(boxes) == 2:
            # Combine geometries if possible when they can be adjusted to an
            # alternative range.
            # Ensure the geometries are ordered from left to right across
            # domain.
            boxes = sorted(boxes, key=lambda box: box.bounds[0])
            invalid_boxes = [box.is_valid for box in boxes]
            invalid_boxes = [
                invalid for invalid in invalid_boxes if not invalid
            ]
            if invalid_boxes:
                print(f"At least one of the two {boxes} is invalid."
                      f"Attempting to fix {boxes}.")
                boxes = [
                    box if box.is_valid else box.buffer(0) for box in boxes
                ]
            if getattr(target_x.units, "modulus", None):
                if utils.ndarray.allclose(
                    boxes[0].bounds[0] % target_x.units.modulus,
                    boxes[1].bounds[2] % target_x.units.modulus,
                ):
                    boxes[0] = shapely.affinity.translate(
                        boxes[0], xoff=target_x.units.modulus
                    )
                    boxes = unary_union(boxes)
```
