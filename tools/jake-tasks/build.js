module.exports = function() {
  var cmds = [
              "uglifyjs --output dist/m.min.js src/m.js"
              ];
  var callback = function() {
  };
  var opts = {
      stdout: true,
      stderr: true,
      breakOnError: false
  };

  jake.exec( cmds, callback, opts );
};
