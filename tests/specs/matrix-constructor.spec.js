var M = require("src/m");

describe("Matrix", function() {
  describe("with dimensions", function() {
    var m;
    beforeEach(function() {
      m = new M.Matrix(2, 2);
    });
    it("should create a new matrix", function() {
      expect(m).toBeDefined();
    });
    it("should initialize new matrix to 0", function() {
      expect(m).toBeLike([0, 0, 0, 0]);
    });
  });
  describe("with initializer", function() {
    var m;
    beforeEach(function() {
      m = new M.Matrix([[1, 2], [3, 4]]);
    });
    it("should create a new matrix", function() {
      expect(m).toBeDefined();
    });
    it("should initialize the new matrix with values", function() {
      expect(m).toBeLike([1, 2, 3, 4]);
    });
  });
});