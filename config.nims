from os import `/`

let root = projectDir()

task test, "Run tests":
  let testDir = root / "test"
  if dirExists(testDir):
    for t in listFiles(testDir):
      echo "Running " & t
      selfExec "c -r --hints:off " & t
