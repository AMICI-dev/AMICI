With `ENABLE_GCOV_COVERAGE`, model sources are no longer compiled with
coverage instrumentation; model extensions are only linked with `--coverage`,
for the instrumented amici library. This allows compiler cache (ccache) hits
for models built in different directories. Coverage of the model sources
themselves is no longer recorded.
