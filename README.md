# opensnoopclose
Opensnoopclose Trace open() and close() syscalls.

Run the `opensnoop_close` with  

`sudo ./opensnoop_close -F`  


and issue a sample dd command for testing the calls  

`dd if=/dev/zero of=./temp_test_file bs=1M count=30 oflag=direct`
