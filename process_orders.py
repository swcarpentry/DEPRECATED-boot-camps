import csv
from itertools import groupby

def main():
    with open('../Downloads/Orders-2953573217.csv', 'r') as f:
        data = [row for row in csv.reader(f)]
        cols, orders = data[0], data[1:]

    index = {col: i for i, col in enumerate(cols)}

    orders = [(int(order[index['Tickets']]), 
               order[index['Date']], 
               order[index['Email Address']]) for order in orders][::-1]

    # Group and sort 
    sorted_orders = []
    for key, grp in groupby(orders, key=lambda x: x[1]):
        sorted_grp = sorted(grp, key=lambda x: x[0], reverse=True) 
        sorted_orders.extend(sorted_grp)
    orders = sorted_orders

    n, tickets = 0, 0
    while tickets < 50:
        tickets += orders[n][0]
        n += 1

    madeit = [order[2] for order in orders[:n]]
    maybe = [order[2] for order in orders[n:n+10]]
    nogo = [order[2] for order in orders[n+10:]]

    print "made it: {}\n".format(', '.join(madeit))
    print "maybe: {}\n".format(', '.join(maybe))
    print "nogo: {}\n".format(', '.join(nogo))


if __name__ == "__main__":
    main()
