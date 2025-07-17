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