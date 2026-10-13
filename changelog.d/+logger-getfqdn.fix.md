AMICI's logger setup no longer calls `socket.getfqdn()`, which could block
for more than 30 seconds per process on hosts whose hostname does not
resolve quickly (e.g., GitHub-hosted macOS 15/26 runners).
