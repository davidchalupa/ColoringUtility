#include <stdio.h>
#include "cli.h"

int main(int argc, char **argv)
{
    cli *cli_instance;

    cli_instance = new cli();
    int cli_return_value = cli_instance->start_cli(argc, argv);
    delete(cli_instance);

    return cli_return_value;
}
