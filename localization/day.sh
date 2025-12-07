#!/bin/bash

# 定义开始和结束日期
start_date="20150101"
end_date="20151231"

# 将开始日期转换为时间戳
start_timestamp=$(date -d "$start_date" +%s)

# 将结束日期转换为时间戳
end_timestamp=$(date -d "$end_date" +%s)

# 逐行生成日期，并写入到文件中
current_timestamp=$start_timestamp
while [ "$current_timestamp" -le "$end_timestamp" ]
do
    # 将当前时间戳转换为日期格式
    current_date=$(date -d "@$current_timestamp" +%Y%m%d)
    
    # 将日期写入到文件中
    echo "$current_date" >> dates.txt

    # 增加一天
    current_timestamp=$(date -d "$current_date + 1 day" +%s)
done

