# Here are my top tips for turning Python into optimized Cython

1. Add @ccall (public fn) and @cfunc (private fn) decorators to each function that should be called directly without runtime module+function lookups.

2. Add type hints to your function inputs and outputs. Remember to use `cy.*` types for C types.

3. Add type hints to all of your local variables. Remember to use `cy.*` types for C types. NOTE: Avoid unboxing values if you are just going to re-box them immediately when they are added to a Python collection.

4. Declare type hints on global variables with `var_name = cy.declare(cy.var_type, var_value)`.

5. Add the `@cy.wraparound(False)` and `@cy.boundscheck(False)` decorators to any functions that perform array access. NOTE: This is not needed for struct/tuple access.

6. Add the `@cy.cdivision(True)` decorator to functions that perform unnecessary divide by zero error checking.

7. Add the `@cy.exceptval(some_c_value)` decorator to any function that might throw an exception. You have to pick `some_c_value` to be something that would never be thrown on success. This will add return value error checking to the call sites where you use this function.

8. When you believe that you can't optimize a function any further, set the `@cy.profile(False)` decorator on it to remove profiling code from it. This will make its calling functions run faster in the profiler and will just merge its runtime into its calling function's runtime.

9. Prefer ctuples (i.e., `tuple[cy.*, cy.*]`) over small Numpy arrays wherever it makes sense. Thus far, this has worked nicely for cell coordinates `coord_yx` and `coord_tyx` and space-time vectors `vec_xy` and `vec_xyz`.

10. Use the `pyidx` type (from `pyretechnics.cy_types` or `pyretechnics.py_types`) for integer indexes into arrays as well as row and column counts.

11. Use `for i in range(n):` to create C `for` loops. Avoid `range(start, stop, step)` calls if possible. Similarly, use `for` loops with accumulator variables instead of Python list comprehensions for speed.

12. Remember to add the function signatures (in Cython syntax) for any public module functions to a file called `mymodule.pxd` when they are implemented in a module called `mymodule.py`. Use `cdef` for C functions and `cpdef` for your `@ccall` decorated functions.

13. Always write conditional imports using `if cython.compiled:` at the top of your module. Provide C functions in the true branch and equivalent Python functions in the false branch.
