from os import `/`
from strutils import endsWith

let root = projectDir()

task test, "Run tests":
  let testDir = root / "test"
  if dirExists(testDir):
    for t in listFiles(testDir):
      if t.endsWith(".nim"):
        echo "Running " & t
        selfExec "c -r --hints:off " & t

# TODO: create documentation
# nim doc --project --index:on --git.url:https://github.com/Psirus/Nimfem --git.commit:master --outdir:htmldocs nimfem.nim
# cp htmldocs/nimfem.html htmldocs/index.html
# ghp-import htmldocs -b gh-pages -p
