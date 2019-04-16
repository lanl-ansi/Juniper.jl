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

To access these statistics you can use `JuMP.backend(m)` i.e.:

```
internal_model = JuMP.backend(m)
println("#IntVars: ", internal_model.optimizer.model.inner.nintvars)
```