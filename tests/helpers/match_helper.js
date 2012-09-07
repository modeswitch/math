function toBeLike(expected) {
  var EPSILON = 0.00001;

  if (typeof(this.actual) == 'number')
    return Math.abs(this.actual - expected) < EPSILON;

  if (this.actual.length != expected.length) return false;
  for (var i = 0; i < this.actual.length; i++) {
    if (isNaN(this.actual[i]) !== isNaN(expected[i]))
      return false;
    if (Math.abs(this.actual[i] - expected[i]) >= EPSILON)
      return false;
  }
  return true;
}

beforeEach(function() {
  this.addMatchers({
    toBeLike: toBeLike
  });
});