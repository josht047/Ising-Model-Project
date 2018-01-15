#!/bin/sh

#this script takes all of the new or modified files in the local repository and commits the changes to the remote git repository on github

#change directory to the local folder
cd /Users/Josh/Desktop/isingmodel 

#switching to add the master branch
git checkout master

#adding all new or modified files
git add -A 

#asking user to input commit message
echo Please type commit message
read message

#committing changes
#commit message goes in commas
git commit -m “message”     

#pushing to remote git repository
git push -u origin master

echo Press Enter to continue
read
