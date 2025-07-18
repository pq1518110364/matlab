# matlab

$ git init    
$ git add .    
$ git commit  

#初始化
git init //创建
git clone /path/to/repository //检出
git config --global user.email "you@example.com" //配置email
git config --global user.name "Name" //配置用户名

#操作
git add <file> // 文件添加，A → B
git add . // 所有文件添加，A → B

git commit -m "代码提交信息" //文件提交，B → C
git commit --amend //与上次commit合并, *B → C

git push origin master //推送至master分支, C → D
git pull //更新本地仓库至最新改动， D → A
git fetch //抓取远程仓库更新， D → C

git log //查看提交记录
git status //查看修改状态
git diff//查看详细修改内容
git show//显示某次提交的内容