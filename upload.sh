version=`date`
echo ${version}
git add -A
git commit -m "${version}"
git branch -M main
git push -u origin main
