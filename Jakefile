var PATH = process.env["PATH"].split( ':' );
var MODULES_BIN = process.cwd() + "/node_modules/.bin";
PATH.unshift( MODULES_BIN );
PATH = PATH.join( ':' );

process.env["PATH"] = PATH;

task( "default", [], require( "./tools/jake-tasks/default" ) );

desc( "lint code" );
task( "lint", [], require( "./tools/jake-tasks/lint" ) );

desc( "compile code and documentation" );
task( "build", ["src/m.js"], require( "./tools/jake-tasks/build" ) );

desc( "remove compiled code" );
task( "clean", [], require( "./tools/jake-tasks/clean" ) );

desc( "run tests" );
task( "test", [], require( "./tools/jake-tasks/test" ) );