import math

from sympy import *

alpha = 0.8
min_best = 1
R = 3
if R < 2:
    print("R should be larger than 2")
    exit(1)
print("alpha: {}".format(alpha))
print("min_best: {}".format(min_best))
print("for yang's R: {}".format(R))
print("from now on, the rank value of each group is reverse order!")
tolerance = 1e-6


# SB scores

def generateStd_SB(group_data):
    y = [float(data[0]) for data in group_data]
    y.sort()
    std_s = min(alpha * y[-1], y[-min_best]) - tolerance
    return std_s, y


def convertRank1(group_data):
    # Start with a score of 0
    cnt = 0
    y = []
    for data in group_data:
        y.append((cnt, float(data[0])))
        cnt += 1
    y.sort(key=lambda x: x[1])
    dict1 = {}
    for i in range(len(y)):
        dict1[y[i][0]] = i
    cnt = 0
    for data in group_data:
        str1 = str(dict1[cnt])
        data[0] = str1
        cnt += 1
    return False


def convertRank2(group_data):
    # 这个是二元分类的mapping 方法
    # 根据y的大小，将y映射到 0 or 1
    std_y, y = generateStd_SB(group_data)
    for data in group_data:
        if float(data[0]) < std_y:
            data[0] = str(0)
        else:
            data[0] = str(1)
    return False


def convertRank3(group_data):
    # yang's method:
    # 根据y的大小，将y映射到 0...R
    std_s, y = generateStd_SB(group_data)
    gap = (std_s - y[0]) / R
    # lst = [(R, std_s, 1e10)]  # rank, lower bound, upper bound
    # for p in range(R - 1, -1, -1):
    #     lst.append((p, std_s - (R - p) * gap, std_s - (R - p - 1) * gap))
    # for data in group_data:
    #     for p in range(R):
    #         if lst[p][1] <= float(data[0]) < lst[p][2]:
    #             data[0] = str(lst[p][0])
    #             break
    # better way to do the above:
    for data in group_data:
        data[0] = str(R - int(max(0., math.ceil((std_s - float(data[0])) / gap))))

    return False


def convertRank4(group_data):
    # 在 std——s 前面的 进行 mapping， 最差的都是0，其余的 按照 比 std_s 大的个数进行 mapping
    std_s, y = generateStd_SB(group_data)
    map1 = {}
    cnt = 0
    for i in range(len(y)):
        if y[i] >= std_s:
            cnt += 1
            map1[y[i]] = cnt

    for data in group_data:
        if float(data[0]) < std_s:
            data[0] = str(0)
        else:
            print(data[0], str(map1[float(data[0])]))
            data[0] = str(map1[float(data[0])])

    return False


def convertRank(group_data, rule):
    if_discard = False
    if rule == '1':
        if_discard = convertRank1(group_data)
    elif rule == '2':
        if_discard = convertRank2(group_data)
    elif rule == '3':
        if_discard = convertRank3(group_data)
    elif rule == '4':
        if_discard = convertRank4(group_data)
    else:
        print("rule : {}. error: no such rule".format(rule))
    return if_discard
