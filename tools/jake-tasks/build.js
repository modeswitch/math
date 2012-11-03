var fs = require("fs");

module.exports = function() {
  var cmds = [
              "cd docs && make html",
              "uglifyjs --output dist/math.min.js src/math.js"
              ];
  var callback = function() {
  };
  var opts = {
      stdout: true,
      stderr: true,
      breakOnError: false
  };

  if(!fs.existsSync("dist")) {
    fs.mkdirSync("dist");
  }
  if(fs.existsSync("dist/math.min.js")) {
    fs.unlinkSync("dist/math.min.js");
  }
  jake.exec( cmds, callback, opts );
};
