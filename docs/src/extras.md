## Statistics

You can get some statistics after the problem is solved.

| Variable  | Description                                  |
|-----------|----------------------------------------------|
| nintvars  | Number of integer variables                  |
| nbinvars  | Number of binary variables                   |
| nnodes    | Number of explored nodes in branch and bound |
| ncuts     | Number of cuts                               |
| nbranches | Number of branches                           |
| nlevels   | Deepest level reached (Root node is level 1) |

To access these statistics you can use `m.internalModel` i.e.:

```
m.internalModel.nintvars
```