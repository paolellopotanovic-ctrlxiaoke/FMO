Orbitals provided here are to be inserted into $FMOHYB,
to be used in (E)FMO runs.
They are provided for one atom, where a fragment boundary is defined.
The orbitals are basis set dependent (but thought to be very weakly affected
by QM method used in making them).
Normally, only single covalent bonds are detachable.

For a X-Y bond there are two asymmetric ways to define the fragment boundary:
(1)
 $FMOBND
     -X     Y 
will define the boundary as X|-Y. You will have to use orbitals for atom X.

On the other hand,
(2)
 $FMOBND
     -Y     X 
will define the boundary as X-|Y. You will have to use orbitals for atom Y.

Which of the two is better depends on the system.

Thus, orbitals provided for C can be used for fragmenting bonds like
C|-C, C|-O, etc.

If you want to make your own orbitals for an atom not provided,
use makeHMO.inp. 
