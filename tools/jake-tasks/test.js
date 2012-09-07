var nodeSpawn = require("child_process").spawn;

spawn = (function(_spawn) { return function(filename, args) {
  // bridge input/output
  // log.info(filename + ' ' + (!!args ? args.join(' ') : ''));
  var proc = _spawn(filename, args);
  process.stdin.pipe(proc.stdin);
  proc.stdout.pipe(process.stdout);
  proc.stderr.pipe(process.stderr);
  proc.on("exit", function(code, signal) {
    process.exit(code);
  });

  return proc;
}; })(nodeSpawn);

module.exports = function() {
  var program = "jasmine-node";
  var args = "--test-dir tests";
  var env = process.env;
  var NODE_PATH = env["NODE_PATH"];
  NODE_PATH = NODE_PATH ? NODE_PATH.split(":") : [];
  NODE_PATH.push(".");
  NODE_PATH = NODE_PATH.join(":");
  env["NODE_PATH"] = NODE_PATH;
  var opts = {
    cwd: process.cwd(),
    env: env
  };

  spawn( program, args.split(" "), opts );
};
