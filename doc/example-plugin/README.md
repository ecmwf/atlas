atlas-example-plugin
====================

This is a Atlas plugin that serves as a template.
The plugin can be used to register concrete implementations for abstract Atlas concepts such as:
Grid, MeshGenerator, Partitioner, FunctionSpace, Interpolation method, ...


Requirements:
-------------
- atlas ... or greater



Loading of atlas plugins:
-------------------------

- If your project explicitly links with the plugin library, then nothing needs to be done.
  The plugin will be loaded as any explicitly linked library.

- If this plugin is installed in the same install-prefix as the executable or the eckit library,
  then nothing needs to be done. The plugin will be automatically detected and loaded at runtime.
  This is the recommended approach.

- For plugins installed in other locations, use the environment variable ATLAS_PLUGIN_PATH

      export ATLAS_PLUGIN_PATH=/path/to/plugin      # colon separated list

  When using older versions of eckit (version <= 1.24.0), instead use

      export PLUGINS_MANIFEST_PATH=/path/to/plugin/share/plugins  # colon separated list
      export LD_LIBRARY_PATH=/path/to/plugin/lib:$LD_LIBRARY_PATH # use DYLD_LIBRARY_PATH for macOS

Verify the plugin is loaded using the `atlas` executable:

    atlas --info

should print a "Plugin" section containing "atlas-example-plugin" with associated version and git sha1.

To debug plugin detection, set environment variable `export MAIN_DEBUG=1`

