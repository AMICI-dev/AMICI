# Changelog fragments

Each pull request with user-visible changes adds its changelog entry as a
separate file to this directory, instead of editing `CHANGELOG.md` directly.
This avoids merge conflicts between pull requests. At release time, the
fragments are collected into `CHANGELOG.md` using
[towncrier](https://towncrier.readthedocs.io/).

## Adding an entry

Create a file `changelog.d/<name>.<type>.md`, where `<type>` is one of

* `breaking`: **BREAKING CHANGES**
* `deprecation`: **Deprecations**
* `feature`: **Features**
* `fix`: **Fixes**
* `performance`: **Performance**
* `doc`: **Documentation**
* `misc`: **Other changes**

and `<name>` is `+` followed by a short descriptive name
(e.g., `+logger-getfqdn.fix.md`). The name only needs to be unique; it does
not appear in `CHANGELOG.md`.

The file contains the text of the entry as Markdown, without a leading
bullet, hard-wrapped at 79 characters. Describe the change from the user's
perspective, and reference relevant issues or pull requests in the text.
For example, `changelog.d/+logger-getfqdn.fix.md`:

```markdown
AMICI's logger setup no longer calls `socket.getfqdn()`, which could block
for more than 30 seconds per process on hosts whose hostname does not
resolve quickly.
```

You can also create the file using `towncrier create`, e.g.,
`towncrier create +logger-getfqdn.fix.md`. towncrier is included in the
`test` dependency group
(`pip install --group python/sdist/pyproject.toml:test`).

To preview the changelog section for the next release, run
`towncrier build --draft --version <next version>` from the repository root.
