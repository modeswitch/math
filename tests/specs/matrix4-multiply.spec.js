var M = require("src/math");

describe("matrix4 multiply", function() {
  describe("with result parameter", function() {
    var m1, m2, result;
    beforeEach(function() {
      m1 = new M.Matrix([[1,2,3,4],[5,6,7,8],[9,1,2,3],[4,5,6,7]]);
      m2 = new M.Matrix([[4,5,6,7],[8,9,1,2],[3,4,5,6],[7,8,9,1]]);
    });
    it("should result the result matrix", function() {
      var result = new M.Matrix(4,4);
      var retval = M.matrix4.multiply(m1, m2, result);
      expect(result).toBeLike([57, 67, 59, 33, 145, 171, 143, 97, 71, 86, 92, 80, 123, 145, 122, 81]);
      expect(retval).toEqual(result);
    });
    it("should update the result matrix", function() {
      var result = new M.Matrix(4,4);
      M.matrix4.multiply(m1, m2, result);
      expect(result).toBeLike([57, 67, 59, 33, 145, 171, 143, 97, 71, 86, 92, 80, 123, 145, 122, 81]);
    });
  });
  describe("without result parameter", function() {
    var m1, m2, result;
    beforeEach(function() {
      m1 = new M.Matrix([[1,2,3,4],[5,6,7,8],[9,1,2,3],[4,5,6,7]]);
      m2 = new M.Matrix([[4,5,6,7],[8,9,1,2],[3,4,5,6],[7,8,9,1]]);
    });
    it("should return the result matrix", function() {
      var retval = M.matrix4.multiply(m1, m2);
      expect(retval).toBeLike([57, 67, 59, 33, 145, 171, 143, 97, 71, 86, 92, 80, 123, 145, 122, 81]);
    });
  });
});