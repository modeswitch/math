var M = require("src/math");

describe("Vector", function() {
  describe("with dimensions", function() {
    var v;
    beforeEach(function() {
      v = new M.Vector(2);
    });
    it("should create a new vector", function() {
      expect(v).toBeDefined();
    });
    it("should initialize new vector to 0", function() {
      expect(v).toBeLike([0, 0]);
    });
  });
  describe("with initializer", function() {
    var v;
    beforeEach(function() {
      v = new M.Vector([1, 2, 3, 4]);
    });
    it("should create a new matrix", function() {
      expect(v).toBeDefined();
    });
    it("should initialize the new matrix with values", function() {
      expect(v).toBeLike([1, 2, 3, 4]);
    });
  });
});