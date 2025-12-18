#!/bin/sh

#echo ">>>>>>>>>>>>>>>> ENTERING COMMAND TEMPLATE <<<<<<<<<<<<<<<<<<<<<<"





echo $$ > .bpipe/commands/1.pid

CMD_SHELL=bash
(command -v bash  > /dev/null) || CMD_SHELL=sh

cat .bpipe/commandtmp/1/cmd_run.sh | setsid $CMD_SHELL -e 

result=$?
echo -n $result > .bpipe/commandtmp/1/cmd.exit.tmp
mv .bpipe/commandtmp/1/cmd.exit.tmp .bpipe/commandtmp/1/cmd.exit



#echo ">>>>>>>>>>>>>>>> EXITING COMMAND TEMPLATE <<<<<<<<<<<<<<<<<<<<<<"

exit $result