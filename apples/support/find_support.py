from apples.support.Bootstrapping import Bootstrapping
import json

def get_support(results):
    support = []
    count = {}
    for result in results:
        query = result[0]["placements"][0]["n"][0]
        count[query] = {}
        count[query]['total'] = 0
        for i in range(1, Bootstrapping.sample_count+1):
            branch = result[i]["placements"][0]["p"][0][0]
            if branch not in count[query]:
                count[query][branch] = 0

            count[query][branch] += 1
            count[query]['total'] += 1

        original_branch = result[0]["placements"][0]["p"][0][0]
        if original_branch not in count[query]:
            support.append(0)
        else:
            support.append(count[query][original_branch]/count[query]['total'])

    # fp = open("placement_count.txt", "w")
    # output = json.dumps(count)
    # fp.write(output)
    # fp.close()

    return support

def get_support_all(results):
    support = {}
    for result in results:
        query = result[0]["placements"][0]["n"][0]
        count = {}
        count['total'] = 0
        for i in range(1, Bootstrapping.sample_count+1):
            branch = result[i]["placements"][0]["p"][0][0]
            if branch not in count:
                count[branch] = 0

            count[branch] += 1
            count['total'] += 1
        
        temp = {}
        for branch in count:
            if branch != 'total':
                temp[branch] = count[branch]/count['total']
                
        support[query] = temp

    return support    