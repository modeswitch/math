var fs = require("fs");

module.exports = function() {
  var cmds = [
              "cd docs && make html",
              "uglifyjs --output dist/m.min.js src/m.js"
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
  if(fs.existsSync("dist/m.min.js")) {
    fs.unlinkSync("dist/m.min.js");
  }
  jake.exec( cmds, callback, opts );
};
