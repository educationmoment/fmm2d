# fmm2d
Fast Marching Method solver

A header-only, first-order Fast Marching Method solver for the Eikonal equation |∇T(x)| F(x) = 1 on uniform 2D grids, following Sethian's Level Set Methods and Fast Marching Methods (2nd ed., Ch. 8).

The solver uses Sethian's upwind finite-difference scheme with Known-only stencil updates and an indexed binary heap with back-pointers, giving O(N log N) complexity. Grids may have unequal spacing (dx ≠ dy), and nodes with F ≤ 0 act as impassable obstacles.
