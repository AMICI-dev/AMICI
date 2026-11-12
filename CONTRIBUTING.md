# Contributing to AMICI

We are happy about contributions to AMICI in any form, be it new functionality,
documentation, bug reports, or anything else.

If you would to contribute to AMICI, a good start is looking for issues tagged
[`good first issue`](https://github.com/AMICI-dev/AMICI/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22)
or [`help wanted`](https://github.com/AMICI-dev/AMICI/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22).
For other ideas or questions, just post an issue.

For code contributions, please read our
[developer's guide](https://amici.readthedocs.io/en/latest/development.html) first.

## Contribution policy

We use AI coding tools ourselves, and we have no objection to contributors
using them.
What we need, with or without AI, is a human who understands the change,
stands behind it, and has put in the work to make it fit the project's
architecture, conventions, and style.
Reviewing a contribution takes real maintainer time —
often more than writing the change would have taken.
A contribution has to be worth more to the project than the time it takes
to review it.
To keep that balance, the following rules apply to all contributions
(code, documentation, issues, and comments).

### Use of AI tools

1. **Disclose any use of AI.**
   State in the pull request or issue description which AI tools you used
   and for what, and which parts of the contribution are AI-generated
   (e.g., "implementation drafted with Claude Code and revised by hand,
   tests written manually").
   Using AI tools to understand the code or to debug does not need to be
   disclosed, nor does routine autocompletion; what matters is what ends up
   in the contribution.
   Pull requests with apparently undisclosed use of AI may be closed
   without review.

2. **You are the author, and you are responsible for the contribution.**
   Before submitting, you should have reviewed the change, understand what
   it does and why, and have run the relevant tests yourself.
   You are also responsible for ensuring that AI-generated content does not
   infringe third-party copyright and can be released under AMICI's license.

3. **We talk to you, not to your agent.**
   You may use AI tools to draft the initial text of an issue or pull
   request, as long as you review it before posting and disclose it.
   In the discussion that follows, however, we expect to talk to you.
   Do not use AI tools to speak for you in review comments or discussions.
   They may help you present information, e.g., format a table or a minimal
   example, but the reasoning and the wording have to be your own.
   If it seems that review comments are answered by relaying them to an AI
   tool and pasting the reply back, we may stop reviewing and close the pull
   request.

### Review capacity

4. **Discuss before you (or your agent) implement.**
   For anything beyond a small, obvious fix (as a rough guide: more than a
   few dozen lines, a new feature, or changes to the code generation or the
   C++ core), open an issue or comment on an existing one first, and wait
   for feedback from a maintainer.
   A small fix can go straight to a pull request; describe the bug there.
   Large pull requests that were not discussed beforehand are closed
   without review.
   This saves time on both sides: yours, if the change is not wanted;
   ours, if it is wanted but needs a different approach.

5. **Contributions have to be worth the review effort.**
   We close issues and pull requests whose benefit to the project does not
   justify the time it takes to review them, regardless of how they were
   produced (bug reports excepted, see below).
   This includes contributions that are mostly AI-generated and do not meet
   the quality standards described in the
   [developer's guide](https://amici.readthedocs.io/en/latest/development.html).

6. **Bug reports are the exception.**
   We would rather have a rough bug report than none.
   Bug reports may be drafted with AI tools, as long as this is disclosed,
   you actually observed the problem, and you tell us how to reproduce it.
   AI-generated analyses of the cause are fine if marked as such, but they
   do not replace a reproducible example: a minimal model or script we can
   run, preferably reduced to the bare minimum rather than your full model
   and data.
   Reports of problems nobody has observed are closed.
   This exception does not extend to feature requests, questions, or
   anything else that starts a discussion.

We are a small team and cannot afford to spend review time on
contributions that mainly serve to pad a GitHub profile.
If we closed your pull request with a link to this section, it is not a
judgement of you or your work; we have to prioritize.
We are glad to support anyone who wants to get involved for longer, brings
some background in the field, and is willing to familiarize themselves
with the code base; we are happy to point you to suitable issues.
If that is you, open an issue and tell us what you would like to work on.
