;(function(window, undefined) {
  'use strict';

  /** Detect free variable `exports` */
  var freeExports = typeof exports == 'object' && exports &&
    (typeof global == 'object' && global && global == global.global && (window = global), exports);

  function nop() {}

  function raise(message) {
    throw new Error(message);
  }

  function isTypedArray(object) {
    return object instanceof ArrayBuffer ||
           object instanceof Int8Array ||
           object instanceof Uint8Array ||
           object instanceof Int16Array ||
           object instanceof Uint16Array ||
           object instanceof Int32Array ||
           object instanceof Uint32Array ||
           object instanceof Float32Array ||
           object instanceof Float64Array;
  }

  function arrayType(object) {
    if(object instanceof ArrayBuffer) {
      return ArrayBuffer;
    } else if(object instanceof Int8Array) {
      return Int8Array;
    } else if(object instanceof Uint8Array) {
      return Uint8Array;
    } else if(object instanceof Int16Array) {
      return Int16Array;
    } else if(object instanceof Uint16Array) {
      return Uint16Array;
    } else if(object instanceof Int32Array) {
      return Int32Array;
    } else if(object instanceof Uint32Array) {
      return Uint32Array;
    } else if(object instanceof Float32Array) {
      return Float32Array;
    } else if(object instanceof Float64Array) {
      return Float64Array;
    } else {
      return undefined;
    }
  }

  function readHeader(a, offset) {
    return a.buffer[offset];
  }

  function writeHeader(a, offset, value) {
    a.buffer[offset] = value;
  }

  // Offsets
  var DIMENSION = 0;
  var TYPE =  1;

  // Types
  var NONE = 0;
  var TRANSFORM = 1;

  // Temporary storage
  var TMP0 = new Float32Array(16);
  var TMP1 = new Float32Array(16);

  /*
   * Constructors
   */
  var ARRAY_TYPE = Float32Array;
  var HEADER_SIZE = 1;  // allocation header size (in elements)
  var ELEMENT_SIZE = ARRAY_TYPE.BYTES_PER_ELEMENT;

  function Matrix(arg0, arg1) {
    var argc = arguments.length;
    var rows, columns, size, values;
    if(1 === argc) {
      // Argument is an array of values
      if(Array.isArray(arg0)) {
        rows = arg0.length;
        columns = arg0[0].length;
        size = rows * columns + HEADER_SIZE;
        values = arg0;
      }
    } else if(2 === argc) {
      // Arguments are rows and columns
      rows = arg0;
      columns = arg1;
      size = rows * columns + HEADER_SIZE;
    } else {
      throw new Error("invalid constructor invocation");
    }

    // NOTE: Allocate an ArrayBuffer to hold a header and data, then create a view
    // for the data only. Internally, we will access the header directly from the buffer,
    // but this way the header stays out of the way for client code.
    var buffer = new ArrayBuffer(size * ELEMENT_SIZE);
    var matrix = new ARRAY_TYPE(buffer, HEADER_SIZE * ELEMENT_SIZE);
    writeHeader(matrix, DIMENSION, columns);
    if(values) {
      for(var i = 0; i < rows; ++ i) {
        var row = values[i];
        if(row.length !== columns) {
          throw new Error("invalid constructor invocation");
        }
        matrix.set(values[i], i * columns);
      }
    }

    return matrix;
  }

  function Transform(dimension) {
    var argc = arguments.length;
    if(1 === argc) {
      // Create an identity matrix large enough to hold an affine transform
      ++ dimension;
      var transform = new Matrix(dimension, dimension);
      writeHeader(transform, TYPE, TRANSFORM);
      for(var i = 0, l = dimension; i < l; ++ i) {
        transform[i * dimension + i] = 1;
      }
      return transform;
    } else {
      throw new Error("invalid constructor invocation");
    }
  }

  function Vector(arg) {
    var argc = arguments.length;
    if(1 === argc) {
      if(Array.isArray(arg) ||
         isTypedArray(arg)) {
        // Argument is an initializer
        var vector = new Matrix(arg.length, 1);
        vector.set(arg);
        return vector;
      } else {
        // Argument is a size
        return new Matrix(arg, 1);
      }
    } else {
      throw new Error("invalid constructor invocation");
    }
  }

  function Quaternion() {
    var argc = arguments.length;
    if(0 === argc) {
      return new Vector(0, 0, 0, 1);
    } else if(4 === argc) {
      // Components
      return new Vector(arguments[0], arguments[1], arguments[2], arguments[3]);
    } else {
      throw new Error("invalid constructor invocation");
    }
  }

  /*
   * Operations
   */

  function add(a1, a2, result) {
    var size = (a1.length === a2.length) ? a1.length : undefined;
    var dimension = (readHeader(a1, DIMENSION) === readHeader(a2, DIMENSION)) ? 
      readHeader(a1, DIMENSION) : undefined;
    if(undefined === size) {
      throw new Error("arrays are not the same size");
    }
    if(undefined === dimension) {
      throw new Error("arrays have mismatched dimensions");
    }
    result = result || new Matrix(size/dimension, dimension);

    // TODO: optimize for size 2, 3, 4, 9, 16

    for(var i = 0, l = size; i < l; ++ i) {
      result[i] = a1[i] + a2[i];
    }

    return result;
  }

  function subtract(a1, a2, result) {
    var size = (a1.length === a2.length) ? a1.length : undefined;
    var dimension = (readHeader(a1, DIMENSION) === readHeader(a2, DIMENSION)) ? 
      readHeader(a1, DIMENSION) : undefined;
    if(undefined === size) {
      throw new Error("arrays are not the same size");
    }
    if(undefined === dimension) {
      throw new Error("arrays have mismatched dimensions");
    }
    result = result || new Matrix(size/dimension, dimension);

    // TODO: optimize for size 2, 3, 4, 9, 16

    for(var i = 0, l = size; i < l; ++ i) {
      result[i] = a1[i] - a2[i];
    }

    return result;
  }

  function equal(a1, a2, e) {
    e = (undefined !== e) ? e : 0.000001;

    // Scalar
    if("number" === typeof a1 && 
       "number" === typeof a2) {
      return Math.abs(a1 - a2) <= e;
    }

    // Check dimensions
    if(readHeader(a1, DIMENSION) !== readHeader(a2, DIMENSION)) {
      return false;
    }

    // TODO: optimize for size 2, 3, 4, 9, 16

    var size = (a1.length === a2.length) ? a1.length : undefined;
    if(undefined === size) {
      return false;
    }
    
    for(var i = 0, l = size; i < l; ++ i) {
      if(Math.abs(a1[i] - a2[i]) > e) {
        return false;
      }
    }

    return true;
  }

  function vector_length(v) {
    var size = v.length;
    var result = 0;

    // TODO: optimize for size 2, 3, 4

    for(var i = 0, l = size; i < l; ++ i) {
      result += v[i] * v[i];
    }

    return Math.sqrt(result);
  }

  function vector_length2(v) {
    var size = v.length;
    var result = 0;

    // TODO: optimize for size 2, 3, 4

    for(var i = 0, l = size; i < l; ++ i) {
      result += v[i] * v[i];
    }

    return result;
  }

  function clear(a, s) {
    s = (undefined === s) ? 0 : s;
    var size = a.length;

    // TODO: optimize for size 2, 3, 4, 9, 16

    for(var i = 0, l = size; i < l; ++ i) {
      a[i] = s;
    }

    return a;
  }

  function scale(a, s, result) {
    var size = a.length;
    var dimension = readHeader(a, DIMENSION);
    result = result || new Matrix(size/dimension, size);

    // TODO: optimize for size 2, 3, 4, 9, 16

    for(var i = 0, l = size; i < l; ++ i) {
      result[i] = result[i] * s;
    }

    return result;
  }

  function vector_uniform(scalar, arg) {
    var size;
    var result;
    if(isTypedArray(arg)) {
      size = arg.length;
      result = arg;
    } else {
      size = arg;
      result = new Vector(size);
    }

    for(var i = 0, l = size; i < l; ++ i) {
      result[i] = scalar;
    }

    return result;
  }

  function matrix_identity(arg) {
    var dimension;
    var result;
    var i, l;
    if(isTypedArray(arg)) {
      dimension = readHeader(arg, DIMENSION);
      if(arg.length/dimension !== dimension) {
        throw new Error("matrix is not square");
      }
      result = arg;
      for(i = 0, l = arg.length; i < l; ++ i) {
        result[i] = 0;
      }
    } else {
      dimension = arg;
      result = new Matrix(dimension, dimension);
    }

    for(var i = 0; i < dimension; ++ i) {
      result[i * dimension + i] = 1;
    }

    return result;
  }

  function quaternion_identity(result) {
    if(result) {
      result[0] = 0;
      result[1] = 0;
      result[2] = 0;
      result[3] = 1;
    } else {
      result = new Quaternion();
    }

    return result;
  }

  function vector_set(v) {
    var size = v.length;
    var argc = arguments.length - 1;

    if(argc < size) {
      throw new Error("insufficient elements for size: " + size);
    }

    // TODO: optimize for size 2, 3, 4

    for(var i = 0, l = size; i < l; ++ i) {
      v[i] = arguments[i+1];
    }

    return v;
  }

  function vector_dot(v1, v2) {
    var size = (v1.length === v2.length) ? v1.length : undefined;
    if(undefined === size) {
      throw new Error("vectors are not the same size");
    }
    var result = 0;

    // TODO: optimize for size 2, 3, 4

    for(var i = 0, l = size; i < l; ++ i) {
      result += v[i] * v[i];
    }

    return result;
  }

  function vector_negate(v, result) {
    var size = v.length;
    result = result || new Vector(size);

    // TODO: optimize for size 2, 3, 4

    for(var i = 0, l = size; i < l; ++ i) {
      result = -1 * v[i];
    }

    return result;
  }

  function vector_limit(v, limit, result) {
    var size = v.length;
    result = result || new Vector(size);

    // TODO: optimize for size 2, 3, 4

    for(var i = 0, l = size; i < l; ++ size) {
      if(v[i] === 0) {
        continue;
      }
      var scale = Math.abs(v[i]);
      var sign = v[i] > 0 ? 1 : -1;
      if(scale > limit) {
        result[i] = sign * limit;
      } else {
        result[i] = v[i];
      }
    }

    return v;
  }

  function clone(a, result) {
    var size = a.length;
    var dimension = readHeader(a, DIMENSION);
    var clone = result || new Matrix(size/dimension, size);

    // TODO: optimize for size 2, 3, 4, 9, 16

    for(var i = 0, l = size; i < l; ++ i) {
      clone[i] = a[i];
    }

    return clone;
  }

  function vector_normalize(v, result) {
    var size = v.length;
    result = result || new Vector(size);
    var length = 0;
    var i;

    // TODO: optimize for size 2, 3, 4

    for(i = 0, l = size; i < l; ++ i) {
      length += v[i] * v[i];
    }

    length = Math.sqrt(length);

    for(i = 0, l = size; i < l; ++ i) {
      result[i] = result[i]/length;
    }

    return result;
  }

  function vector_distance(v1, v2) {
    var size = (v1.length === v2.length) ? v1.length : undefined;
    if(undefined === size) {
      throw new Error("vectors are not the same size");
    }

    // TODO: optimize for size 2, 3, 4

    var d;
    var r = 0;
    for(var i = 0, l = size; i < l; ++ i) {
      d = v1[i] - v2[i];
      r += d * d;
    }

    return Math.sqrt(r);
  }

  function vector_angle(v1, v2) {
    var size = (v1.length === v2.length) ? v1.length : undefined;
    if(undefined === size) {
      throw new Error("vectors are not the same size");
    }
    var result = 0;

    // TODO: optimize for size 2, 3, 4

    for(var i = 0, l = size; i < l; ++ i) {
      result += v1[i] * v2[i];
    }

    return Math.acos(result);
  }

  function vector_lerp(v1, v2, s, result) {
    var size = (v1.length === v2.length) ? v1.length : undefined;
    if(undefined === size) {
      throw new Error("vectors are not the same size");
    }
    result = result || new Vector(size);

    // TODO: optimize for size 2, 3, 4

    for(var i = 0, l = size; i < l; ++ i) {
      result[i] = v1[i] + s * (v2[i] - v1[i]);
    }

    return result;
  }

  function vector_extract(v, offset, length, result) {
    result = result || new Vector(length);
    for(var i = offset, l = length; i < l; ++ i) {
      result[i] = v[i];
    }
  }

  function vector_direction(v1, v2, result) {
    var size = (v1.length === v2.length) ? v1.length : undefined;
    if(undefined === size) {
      throw new Error("vectors are not the same size");
    }
    result = result || new Vector(size);
    var length = 0;

    // TODO: optimize for size 2, 3, 4

    var i;
    for(i = 0, l = size; i < l; ++ i) {
      result[i] = v2[i] - v1[i];
      length += result[i] * result[i];
    }

    length = 1/Math.sqrt(length);

    for(i = 0, l = size; i < l; ++ i) {
      result[i] = result[i] * length;
    }

    return result;
  }

  function scalar_clamp(s, min, max) {
    s = (s < min) ? min : s;
    s = (s > max) ? max : s;
    return s;
  }

  function vector3_cross(v1, v2, result) {
    result = result || new Vector(3);

    var v1_0 = v1[0],
        v1_1 = v1[1],
        v1_2 = v1[2];
    var v2_0 = v2[0],
        v2_1 = v2[1],
        v2_2 = v2[2];

    result[0] = (v1_1 * v2_2) - (v2_1 * v1_2);
    result[1] = (v1_2 * v2_0) - (v2_2 * v1_0);
    result[2] = (v1_0 * v2_1) - (v2_0 * v1_1);

    return result;
  }

  function vector3_unproject(v, view, projection, viewport, result) {

  }

  function toMathML(a) {
    var size = a.length;   
    var columns = readHeader(a, DIMENSION);
    var rows = size/columns;
    var result = "<mfenced><mtable>";
    for(var i = 0; i < rows; ++ i) {
      result += "<mtr>";
      for(var j = 0; j < columns; ++ j) {
        result += "<mtd><mn>";
        result += a[i * columns + j];
        result += "</mn></mtd>";
      }
      result += "</mtr>";
    }
    result += "</mtable></mfenced>";

    return result;
  }

  function toString(a) {

  }

  function vector_unit(offset, size, result) {
    if(!result) {
      result = new Vector(size);
    } else if(result.length !== size) {
      throw new Error("result has incorrect length " + result.length + ", expected " + size);
    }

    // TODO: optimize for size 2, 3, 4

    for(var i = 0, l = size; i < l; ++ i) {
      result[i] = 0;
    }
    result[offset] = 1;

    return result;
  }

  function quaternion_inverse(q, result) {
    result = result || new Quaternion();

    var q0 = q[0], q1 = q[1], q2 = q[2], q3 = q[3];
    var dot = q0*q0 + q1*q1 + q2*q2 + q3*q3;
    var idot;

    if(0 === dot) {
      result[0] = 0;
      result[1] = 0;
      result[2] = 0;
      result[3] = 0;
    } else {
      idot = 1/dot;
      result[0] = -q0*idot;
      result[1] = -q1*idot;
      result[2] = -q2*idot;
      result[3] = -q3*idot;
    }

    return result;
  }

  function quaternion_conjugate(q, result) {
    result = result || new Quaternion();

    result[0] = -q[0];
    result[1] = -q[1];
    result[2] = -q[2];
    result[3] = q[3];

    return result;
  }

  function matrix_transpose(m, result) {
    var size = m.length;
    var dimension = readHeader(m, DIMENSION);
    var tmp;

    if(size/dimension !== dimension) {
      throw new Error("matrix is not square");
    }

    // TODO: implement algorithm for non-square matrices

    // FIXME: this is not a high-performance algorithm
    for(var i = 0, l = dimension - 2; i <= l; ++ i) {
      for(var j = i+1, k = dimension - 1; j <= k; ++ j) {
        tmp = m[i * dimension + j];
        m[i * dimension + j] = m[j * dimension + i];
        m[j * dimension + i] = tmp;
      }
    }
  }

  function matrix_multiply(m1, m2, result) {
    var m1_size = m1.length;
    var m2_size = m2.length;
    var m1_dimension = readHeader(m1, DIMENSION);
    var m2_dimension = readHeader(m2, DIMENSION);

    if(m1_dimension !== m2_size/m2_dimension) {
      throw new Error("arguments have mismatched rows and columns");
    }

    // TODO: optimize for size 4, 9, 16

    var m = m1_size/m1_dimension;
    var p = m1_dimension;
    var n = m2_dimension;
    var left, right;

    if(result) {
      var result_size = result.length;
      var result_dimension = readHeader(result, DIMENSION);
      if(result_size !== m * n ||
         result_dimension !== n) {
        throw new Error("result has mismatched dimensions");
      }
      if(m1 === result) {
        TMP0.set(m1);
        left = TMP0;
        right = m2;
      } else if(m2 === result) {
        left = m1;
        TMP0.set(m2);
        right = TMP0;
      }
    } else {
      result = new Matrix(m, n);
      left = m1;
      right = m2;
    }

    var tmp;
    for(var i = 0; i < m; ++ i) {
      for(var j = 0; j < n; ++ j) {
        tmp = 0;
        for(var k = 0; k < p; ++k) {
          tmp += left[p * i + k] * right[n * k + j];
        }
        result[n * i + j] = tmp;
      }
    }

    return result;
  }

  function matrix_inverse(m, result) {
    var size = m.length;
    var dimension = readHeader(m, DIMENSION);

    if(size/dimension !== dimension) {
      throw new Error("matrix is not square");
    }

    // TODO: implement inverse for NxN matrices

    result = result || new Matrix(dimension, dimension);
    var a00, a01, a02, a03,
        a10, a11, a12, a13,
        a20, a21, a22, a23,
        a30, a31, a32, a33,

        b00, b01, b02, b03,
        b04, b05, b06, b07,
        b08, b09, b10, b11;
    var determinant, inverseDeterminant;

    if(3 === dimension) {
      a00 = m[0]; a01 = m[1]; a02 = m[2];
      a10 = m[3]; a11 = m[4]; a12 = m[5];
      a20 = m[6]; a21 = m[7]; a22 = m[8];

      b00 = a22 * a11 - a12 * a21;
      b01 = -a22 * a10 + a12 * a20;
      b02 = a21 * a10 - a11 * a20;

      determinant = a00 * b00 + a01 * b01 + a02 * b02;
      if(!determinant) {
        throw new Error("matrix is not invertible");
      }
      inverseDeterminant = 1/determinant;

      result[0] = b00 * inverseDeterminant;
      result[1] = (-a22 * a01 + a02 * a21) * inverseDeterminant;
      result[2] = (a12 * a01 - a02 * a11) * inverseDeterminant;
      result[3] = b01 * inverseDeterminant;
      result[4] = (a22 * a00 - a02 * a20) * inverseDeterminant;
      result[5] = (-a12 * a00 + a02 * a10) * inverseDeterminant;
      result[6] = b02 * inverseDeterminant;
      result[7] = (-a21 * a00 + a01 * a20) * inverseDeterminant;
      result[8] = (a11 * a00 - a01 * a10) * inverseDeterminant;
    } else if(4 === dimension) {
      a00 = m[0]; a01 = m[1]; a02 = m[2]; a03 = m[3];
      a10 = m[4]; a11 = m[5]; a12 = m[6]; a13 = m[7];
      a20 = m[8]; a21 = m[9]; a22 = m[10]; a23 = m[11];
      a30 = m[12]; a31 = m[13]; a32 = m[14]; a33 = m[15];

      b00 = a00 * a11 - a01 * a10;
      b01 = a00 * a12 - a02 * a10;
      b02 = a00 * a13 - a03 * a10;
      b03 = a01 * a12 - a02 * a11;
      b04 = a01 * a13 - a03 * a11;
      b05 = a02 * a13 - a03 * a12;
      b06 = a20 * a31 - a21 * a30;
      b07 = a20 * a32 - a22 * a30;
      b08 = a20 * a33 - a23 * a30;
      b09 = a21 * a32 - a22 * a31;
      b10 = a21 * a33 - a23 * a31;
      b11 = a22 * a33 - a23 * a32;

      determinant = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
      if(!determinant) {
        throw new Error("matrix is not invertible");
      }
      inverseDeterminant = 1/determinant;

      result[0] = (a11 * b11 - a12 * b10 + a13 * b09) * inverseDeterminant;
      result[1] = (-a01 * b11 + a02 * b10 - a03 * b09) * inverseDeterminant;
      result[2] = (a31 * b05 - a32 * b04 + a33 * b03) * inverseDeterminant;
      result[3] = (-a21 * b05 + a22 * b04 - a23 * b03) * inverseDeterminant;
      result[4] = (-a10 * b11 + a12 * b08 - a13 * b07) * inverseDeterminant;
      result[5] = (a00 * b11 - a02 * b08 + a03 * b07) * inverseDeterminant;
      result[6] = (-a30 * b05 + a32 * b02 - a33 * b01) * inverseDeterminant;
      result[7] = (a20 * b05 - a22 * b02 + a23 * b01) * inverseDeterminant;
      result[8] = (a10 * b10 - a11 * b08 + a13 * b06) * inverseDeterminant;
      result[9] = (-a00 * b10 + a01 * b08 - a03 * b06) * inverseDeterminant;
      result[10] = (a30 * b04 - a31 * b02 + a33 * b00) * inverseDeterminant;
      result[11] = (-a20 * b04 + a21 * b02 - a23 * b00) * inverseDeterminant;
      result[12] = (-a10 * b09 + a11 * b07 - a12 * b06) * inverseDeterminant;
      result[13] = (a00 * b09 - a01 * b07 + a02 * b06) * inverseDeterminant;
      result[14] = (-a30 * b03 + a31 * b01 - a32 * b00) * inverseDeterminant;
      result[15] = (a20 * b03 - a21 * b01 + a22 * b00) * inverseDeterminant;
    } else {
      throw new Error("matrix must have dimension 3 or 4");
    }

    return result;
  }



  var api = {
    // Configuration
    debug: false,

    // Utility functions
    isTypedArray: isTypedArray,
    arrayType: arrayType,

    // Constructors
    Matrix: Matrix,
    Transform: Transform,
    Vector: Vector,
    Quaternion: Quaternion,
    Frustum: undefined,
    Perspective: undefined,
    Orthographic: undefined,
    LookAt: undefined,      

    // Constants
    PI: Math.PI,   
    TAU: Math.PI * 2,
    
    // Operators
    clone: clone,
    add: add,
    subtract: subtract,
    equal: equal,
    clear: clear,
    scale: scale,
    toMathML: toMathML,
    toString: toString,
    scalar: {
      clamp: scalar_clamp
    },
    vector: {
      length: vector_length,
      length2: vector_length2,
      set: vector_set,
      dot: vector_dot,
      negate: vector_negate,
      limit: vector_limit,
      normalize: vector_normalize,
      distance: vector_distance,
      angle: vector_angle,
      lerp: vector_lerp,
      direction: vector_direction,
      extract: vector_extract,
      zero: vector_uniform.bind(undefined, 0),
      one: vector_uniform.bind(undefined, 1)
    },
    vector2: {
      x: vector_unit.bind(undefined, 0, 2),
      y: vector_unit.bind(undefined, 1, 2),
      u: vector_unit.bind(undefined, 0, 2),
      v: vector_unit.bind(undefined, 1, 2)
    },
    vector3: {
      cross: vector3_cross,
      unproject: undefined,
      x: vector_unit.bind(undefined, 0, 3),
      y: vector_unit.bind(undefined, 1, 3),
      z: vector_unit.bind(undefined, 2, 3)
    },
    vector4: {
      x: vector_unit.bind(undefined, 0, 4),
      y: vector_unit.bind(undefined, 1, 4),
      z: vector_unit.bind(undefined, 2, 4),
      w: vector_unit.bind(undefined, 3, 4)
    },
    quaternion: {
      identity: quaternion_identity,
      inverse: quaternion_inverse,
      conjugate: quaternion_conjugate,
      slerp: undefined,
      // Quaternion of rotation between two vectors
      rotation: undefined,
      // Convert between quaternion and axis-angle
      toAxisAngle: undefined,
      fromAxisAngle: undefined
    },
    matrix: {      
      multiply: matrix_multiply,
      inverse: matrix_inverse,
      transpose: matrix_transpose,  // Square matrix only!
      set: undefined,
      get: undefined,
      identity: matrix_identity
    },
    matrix2: {
      // Convert between rotation matrix and angle
      toAngle: undefined,
      fromAngle: undefined
    },
    matrix3: {
      // Convert between rotation matrix and quaternion
      toQuaternion: undefined,
      fromQuaternion: undefined,
      // Convert between rotation matrix and axis-angle
      toAxisAngle: undefined,
      fromAxisAngle: undefined
    },
    transform: {
      translate: undefined,
      rotate: undefined,
      scale: undefined,
      // Singular value decomposition
      svd: undefined,
      // Extract the linear part of the transform
      linear: undefined
    }
  };

/*--------------------------------------------------------------------------*/

  // expose the API (borrowed from Lodash: https://raw.github.com/bestiejs/lodash/v0.6.1/lodash.js)
  // some AMD build optimizers, like r.js, check for specific condition patterns like the following:
  if (typeof define == 'function' && typeof define.amd == 'object' && define.amd) {
    // Expose api to the global object even when an AMD loader is present in
    // case Lo-Dash was injected by a third-party script and not intended to be
    // loaded as a module. The global assignment can be reverted in the Lo-Dash
    // module via its `noConflict()` method.
    window.M = api;

    // define as an anonymous module so, through path mapping, it can be
    // referenced as the "underscore" module
    define(function() {
      return api;
    });
  }
  // check for `exports` after `define` in case a build optimizer adds an `exports` object
  else if (freeExports) {
    // in Node.js or RingoJS v0.8.0+
    if (typeof module == 'object' && module && module.exports == freeExports) {
      (module.exports = api).M = api;
    }
    // in Narwhal or RingoJS v0.7.0-
    else {
      freeExports.M = api;
    }
  }
  else {
    // in a browser or Rhino
    window.M = api;
  }
}(this));