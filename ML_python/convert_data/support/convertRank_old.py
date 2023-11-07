from sympy import *

alpha = 0.8
# min_one = 2
# alpha = 0.5
print("alpha: {}".format(alpha))
min_one = 2
# if the ratio of num between s^* and std_s is lager than max_one_alpha * len(group), we discard this group
max_one_alpha = 0.5
# max_one_alpha = 0.9
print("max_one_alpha: {}".format(max_one_alpha))

Tolerance = 1e-6

# for yang's method:
# R = 20
R = 3
epsilon = 1e-15

# for my method:
beta_tolerance = 1e-3
my_epsilon = 1e-4


# for my second method:


def generateStd_SB(group_data, tolerance):
    y = [float(data[0]) for data in group_data]
    y.sort(reverse = True)
    std_s = min(alpha * y[0], y[min_one - 1]) - tolerance
    max_one = int(max_one_alpha * len(group_data))
    if_discard = False
    if y[max_one] > std_s:
        if_discard = True
    return std_s, y, if_discard


def convertRank1(group_data):
    # 这个是递增排序的mapping 方法
    # 根据y的大小，将y映射到1....
    cnt = 0
    y = []
    for data in group_data:
        y.append((cnt, float(data[0])))
        cnt += 1
    y.sort(key = lambda x: x[1], reverse = True)
    dict1 = {}
    for i in range(len(y)):
        dict1[y[i][0]] = i + 1
    cnt = 0
    for data in group_data:
        str1 = str(dict1[cnt])
        data[0] = str1
        cnt += 1
    return False


def convertRank2(group_data):
    # 这个是二元分类的mapping 方法
    # 根据y的大小，将y映射到 1 or 2
    std_y, y, if_discard = generateStd_SB(group_data, Tolerance)
    if if_discard:
        return True
    for data in group_data:
        if float(data[0]) < std_y:
            data[0] = str(2)
        else:
            data[0] = str(1)
    return False


def convertRank3(group_data):
    # yang's method:
    # 根据y的大小，将y映射到 1...R
    std_s, y, if_discard = generateStd_SB(group_data, epsilon)
    if if_discard:
        return True
    gap = (std_s - y[- 1]) / (R - 1) + epsilon
    lst = [(1, std_s, 1e10)]  # rank, lower bound, upper bound
    for p in range(1, R):
        lst.append((p + 1, std_s - p * gap, std_s - (p - 1) * gap))
    for data in group_data:
        for p in range(R):
            if lst[p][1] - Tolerance < float(data[0]) <= lst[p][2]:
                data[0] = str(lst[p][0])
                break
    return False


def convertRank4(group_data):
    std_s, y, if_discard = generateStd_SB(group_data, my_epsilon)
    if if_discard:
        return True
    k = 1 / (log(std_s + beta_tolerance - my_epsilon) - log(y[0] + beta_tolerance + my_epsilon))
    a = 1 - k * log(y[0] + beta_tolerance + my_epsilon)
    # count how many 1s in the group
    # cnt = 0
    # ls = []
    for data in group_data:
        # old_data = data[0]
        data[0] = str(max((floor(k * log(float(data[0]) + beta_tolerance) + a)), floor(1)))
        # if data[0] == '1':
        #     cnt += 1
        #     ls.append(old_data)
    # pprint("cnt: {}".format(cnt))
    # if cnt >= 4:
    #     pprint("ls: {}".format(ls))
    return False


def convertRank5(group_data):
    k = 5
    std_s, y, if_discard = generateStd_SB(group_data, my_epsilon)
    if if_discard:
        return True
    best_score = y[0] + beta_tolerance + my_epsilon
    second_score = std_s - my_epsilon
    # count how many candidates between best_score and second_score
    map1 = {}
    cnt = 1
    for i in y:
        if second_score <= i <= best_score:
            map1[i] = str(cnt)
            cnt += 1
    cnt = max(cnt, 100)
    b = (cnt - 1) / (best_score ** k - second_score ** k)
    a = 1 + b * best_score ** k
    for data in group_data:
        if float(data[0]) >= second_score:
            data[0] = map1[float(data[0])]
            print("data: {}".format(data))
        else:
            rank = (-b * float(data[0]) ** k + a)
            data[0] = str(max(1, int(floor(rank))))

    map2 = {}
    for data in group_data:
        map2[data[0]] = 0
    for data in group_data:
        map2[data[0]] += 1
    # sort map by its key
    map2 = sorted(map2.items(), key = lambda x: int(x[0]))
    pprint("map2: {}".format(map2))
    exit(0)
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
    elif rule == '5':
        if_discard = convertRank5(group_data)
    else:
        print("rule : {}. error: no such rule".format(rule))
    return if_discard
