{
	if (line == 0) {
		for (i = 1; i <= NF; ++i) {
			header[i] = $i;
		}
		headers = NF;
	} else {
		if ($1 == "Donor")
			file = prefix ".Donor.tab";
		else
			file = prefix "." substr($1, 1, 2) ".tab";
		if (first[file] == 0) {
			first[file] = 1;
			printf("%s", header[1]) > file;
			for (i = 2; i <= headers; ++i)
				printf("\t%s", header[i]) >> file;
			printf("\n") >> file;
		}
		print $0 >> file;
	}
	line++;
}
