# wmaee

A new generation of `wmaee`, version 2.0...

## Basics of git
### Merging
To merge branch `david` into `wmaee2`:
1. `git checkout wmaee2`: switch to branch `wmaee2`
2. `git merge david`: merge david into `wmaee2`
3. `git push`: sync local changes with the server
4. `git checkout david`: switch back to `david` for further developments
5. `git rebase wmaee2`: sync branch `david` with the latest `wmaee2` content

`git push --set-upstream origin wmaee2
