eval "$(ssh-agent -s)"
echo "Deploying to marketplace"
chmod 600 .travis/travis_rsa
ssh-add .travis/travis_rsa
git remote add deploy git@danse.chem.utk.edu:/home/git/sasmodels
git push deploy master
