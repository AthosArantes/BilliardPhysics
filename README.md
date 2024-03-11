# Billiard Physics
A library for calculating billiard physics. This project takes the physics code from FooBillard++ into a modern design.

**Changes**

 - Renamed classes/struct for better clarification;
 - Colliders are groups of shapes, each shape/collider have it's own AABB for better performance;
 - Added pocket edge collisions (Prevents ball from clipping through the table when falling into the pocket at very slow speeds);
 - The table "slate" or "cloth" is the plane zero, not a collider;

**Requirements**

C++ 17 Compiler and the std library (vector and algorithm).

**Notes**

The interpretation of moved balls, collisions (ball-cushion, ball-ball) is still up to the developer to implement.

## FooBillard++ License

Copyright (C) 2001 Florian Berger (foobillard)

Copyright (C) 2010/2011 Holger Schaekel (foobillard++)

email:  [foobillardplus@go4more.de](mailto:foobillardplus@go4more.de)

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License Version 2 as published by the Free Software Foundation;

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA