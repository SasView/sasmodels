if [ "$encrypted_cb04388797b6_iv" ]
then
    eval "$(ssh-agent -s)"
    chmod 600 .travis/travis_rsa
    ssh-add .travis/travis_rsa
    git remote add deploy git@danse.chem.utk.edu:/home/git/sasmodels
    # Travis is reporting: [remote rejected] master -> master (shallow update not allowed)
    # Maybe the following will fix it?
    # https://stackoverflow.com/questions/28983842/remote-rejected-shallow-update-not-allowed-after-changing-git-remote-url#28985327
    git fetch --unshallow deploy
    git push deploy master
fi