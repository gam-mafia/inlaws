# Starting unrelated work

Unless the user instructs otherwise, start new work that is unrelated to the
current task in a clean worktree or on a clean branch based on `main`, even when
the request arrives in the same chat. Keep existing unrelated changes separate
rather than carrying them into the new task.

If it is unclear whether a request belongs to the current task, or which branch
or worktree should be used, ask the user to confirm how to proceed before
starting the new work.

# Pull requests

Open a pull request only in response to an explicit request from the user.
Finishing development work, committing changes, merging, or pushing does not
by itself authorize creating a pull request.

# Maintaining the pkgdown reference index

When adding documentation topics, add them to the manually curated reference
index under `reference:` in `_pkgdown.yml`. Keep the index up to date when topics
are renamed or removed. Topics deliberately marked as internal do not need an
index entry.
