from os import `/`

let root = projectDir()

task test, "Run tests":
  let testDir = root / "test"
  if dirExists(testDir):
    for t in listFiles(testDir):
      selfExec "c -r " & t
