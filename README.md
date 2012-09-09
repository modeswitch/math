# MathJS

[![Build Status](https://secure.travis-ci.org/gladiusjs/mathjs.png?branch=develop)](http://travis-ci.org/gladiusjs/mathjs)

MathJS is a math and geometry library for JavaScript. It's optimized for use in high-performance computing applications like 3D games, and it can handle general-purpose operations on matrices and vectors.

## Math for Games

Common algorithms for 2D and 3D games are supported, including affine transform matrices, transforming points or directions, and optimized multiplication for square matrices in 2, 3, and 4 dimensions.

````javascript
// Create a new 3D transform (a 4x4 matrix)
var t = new M.Transform(3);

// Translate and store the result back in t
M.transform.translate(t, [1, 2, 0], t);

// Pre-translate by by changing the parameter order
M.transform.translate([3, 0, 2], t, t);

// Create a new 3D vector
var p = new M.Vector([1, 1, 1]);

// Transform p
M.vector3.transformPoint(t, p, p);
````

## Examples

1. Create an 8x8 identity matrix

        var i = M.matrix.identity(8);

    or

        var i = new M.Matrix(8, 8);
        M.matrix.identity(i);

2. Multiply two 3x3 matrices

        var a = new M.Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
        var b = new M.Matrix([[9, 8, 7], [6, 5, 4], [3, 2, 1]]);
        var product = M.matrix3.multiply(a, b);

3. Use optimized algorithms to multiply two transforms

        var t1 = new M.Transform(3)
        var t2 = new M.Transform(3);
        M.matrix4.multiply(t1, t2, t1);

4. Create a 2D rotation matrix

        var TMP = new Matrix(2, 2);
        var rotation = M.matrix2.fromAngle(M.TAU/2, TMP);

## Reference

Documentation for all functionality and conventions is available [here](http://mathjs.readthedocs.org).