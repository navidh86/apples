from apples.support.Bootstrapping import Bootstrapping


def get_support(results):
    support = []
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

        original_branch = result[0]["placements"][0]["p"][0][0]
        if original_branch not in count:
            support.append(0)
        else:
            support.append(count[original_branch]/count['total'])

    return support