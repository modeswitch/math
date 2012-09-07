module.exports = function() {
  var cmds = [
              "NODE_PATH=$NODE_PATH:. jasmine-node --test-dir tests"
              ];
  var callback = function() {
  };
  var opts = {
      printStdout: true,
      printStderr: true,
      breakOnError: false
  };

  jake.exec( cmds, callback, opts );
};
