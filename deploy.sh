openssl aes-256-cbc -K $encrypted_fe6026add10a_key -iv $encrypted_fe6026add10a_iv -in .travis/travis_rsa.enc -out .travis/travis_rsa
eval "$(ssh-agent -s)"
chmod 600 .travis/travis_rsa
ssh-add .travis/travis_rsa
git remote add deploy git@danse.chem.utk.edu:/home/git/sasmodels
git push deploy master