if( typeof define !== "function" ) {
  var define = require( "amdefine" )( module );
}

define(function() {
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

  /*
   * Constructors
   */
  function Matrix() {
    var argc = arguments.length;
    if(2 === argc) {
      // If we have two values, it must be rows and columns
      return new Float32Array(arguments[0] * arguments[1]);
    } else {
      throw new Error("invalid constructor invocation");
    }
  }

  function Transform(dimensions) {
    // Create an identity matrix large enough to hold an affine transform
    ++ dimensions;
    var transform = new Matrix(dimensions, dimensions);
    for(var i = 0, l = dimensions; i < l; ++ i) {
      transform[i * dimensions + i] = 1;
    }
    return transform;
  }

  function Vector() {
    var argc = arguments.length;
    if(1 === argc) {
      // If a single value is passed, it must be the size of the vector
      return new Float32Array(arguments[0]);
    } else if(argc > 1) {
      // If we received mutliple values, they must be the initializers for the vector
      var vector = new Float32Array(argc);
      for(var i = 0, l = argc; i < l; ++ i) {
        vector[i] = arguments[i];
      }
      return vector;
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
    if(undefined === size) {
      throw new Error("vectors are not the same size");
    }
    result = result || new Float32Array(size);

    // TODO: optimize for size 2, 3, 4, 9, 16

    for(var i = 0, l = size; i < l; ++ i) {
      result[i] = a1[i] + a2[i];
    }

    return result;
  }

  function subtract(a1, a2, result) {
    var size = (a1.length === a2.length) ? a1.length : undefined;
    if(undefined === size) {
      throw new Error("vectors are not the same size");
    }
    result = result || new Float32Array(size);

    // TODO: optimize for size 2, 3, 4, 9, 16

    for(var i = 0, l = size; i < l; ++ i) {
      result[i] = a1[i] - a2[i];
    }

    return result;
  }

  function equal(a1, a2, e) {
    e = (undefined !== e) ? e : 0.000001;
    if(!isNaN(a1) && !isNaN(a2)) {
      return Math.abs(a1 - a2) <= e;
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
    result = result || new Float32Array(size);

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

  function vector_set(v) {
    var size = v.length;
    var argc = arguments.length - 1;

    if(argc < size) {
      throw new Error("insufficient elements for size: " + size);
    }

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

    for(var i = 0, l = size; i < l; ++ i) {
      result += v[i] * v[i];
    }

    return result;
  }

  function vector_negate(v, result) {
    var size = v.length;
    result = result || new Vector(size);

    for(var i = 0, l = size; i < l; ++ i) {
      result = -1 * v[i];
    }

    return result;
  }

  function vector_limit(v, limit, result) {
    var size = v.length;
    result = result || new Vector(size);

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
    var clone = result || new Float32Array(size);

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

  function scalar_equal(s1, s2, e) {
    return Math.abs(s1 - s2) > e;
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

  return {
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
    π: Math.PI,
    TAU: Math.PI * 2,
    τ: Math.PI * 2,
    
    // Operators
    clone: clone,
    add: add,
    subtract: subtract,
    equal: equal,
    clear: clear,
    scale: scale,
    scalar: {
      clamp: scalar_clamp,
      equal: scalar_equal
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
      x: undefined,
      y: undefined,
      u: undefined,
      v: undefined
    },
    vector3: {
      cross: vector3_cross,
      unproject: undefined,
      x: undefined,
      y: undefined,
      z: undefined
    },
    vector4: {
      x: undefined,
      y: undefined,
      z: undefined,
      w: undefined
    },
    quaternion: {
      identity: undefined,
      inverse: undefined,
      conjugate: undefined,
      slerp: undefined,
      // Quaternion of rotation between two vectors
      rotation: undefined, 
      // Convert between quaternion and axis-angle
      toAxisAngle: undefined,
      fromAxisAngle: undefined
    },
    matrix: {
      multiply: undefined,
      inverse: undefined,
      transpose: undefined,
      determinant: undefined,
      get: undefined,
      set: undefined,
      extract: undefined,
      identity: undefined
    },
    matrix2: {
      multiply: undefined,
      // Convert between rotation matrix and angle
      toAngle: undefined,
      fromAngle: undefined
    },
    matrix3: {
      multiply: undefined,
      // Convert between rotation matrix and quaternion
      toQuaternion: undefined,
      fromQuaternion: undefined,
      // Convert between rotation matrix and axis-angle
      toAxisAngle: undefined,
      fromAxisAngle: undefined
    },
    matrix4: {
      multiply: undefined,
    },
    transform: {
      translate: undefined,
      pretranslate: undefined,
      rotate: undefined,
      prerotate: undefined,
      scale: undefined,
      prescale: undefined,
      // Singular value decomposition
      svd: undefined,
      // Extract the linear part of the transform
      linear: undefined
    }
  };
});