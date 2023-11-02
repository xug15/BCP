version=`date`
echo ${version}
git add *
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:xug15/BCP.git
#git push -u origin main
