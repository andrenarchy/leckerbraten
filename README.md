# leckerbraten


Leckerbraten is a tool for simulating heat distribution in solid objects.

Leckerbraten is a tool for simulating heat distribution in solid objects. The
initial goal was to simulate a roast (in german: "Braten") in the oven.
"Lecker" means delicious. :)

If you are interested in the internals of leckerbraten, here is how it
works. Leckerbraten approximates the solution of the [heat
equation](http://en.wikipedia.org/wiki/Heat_equation) (with non-constant
thermal diffusivity) in 2D or 3D with the [finite element method
(FEM)](http://en.wikipedia.org/wiki/Finite_element_method). The code is
based on [http://fenicsproject.org/](FEniCS/Dolfin).

## Dependencies
* [dolfin](http://fenicsproject.org/)
