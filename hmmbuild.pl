open(ALN, "ls /home/dustin/projects/ian_chat/smart_01_06_2016/aln |");

while (<ALN>) {
    chomp;
    $name = $_;
    $name =~ s/aln/hmm/g;
    $command = "hmmbuild $name /home/dustin/projects/ian_chat/smart_01_06_2016/aln/$_";
    system("$command");
}

close(ALN);
