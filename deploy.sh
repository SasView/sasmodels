if [ "$encrypted_cb04388797b6_iv" ]
then
    eval "$(ssh-agent -s)"
    chmod 600 .travis/travis_rsa
    ssh-add .travis/travis_rsa
    git remote add deploy git@danse.chem.utk.edu:/home/git/sasmodels
    git push deploy master
fi
