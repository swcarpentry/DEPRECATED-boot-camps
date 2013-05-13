BEGIN { total = 0 }
/ATOM/ {
    print $3, $5, $6, $7;
    total += 1
}
/AUTHOR/ { print $2, $3 }
END { print "Total atoms: " total }
