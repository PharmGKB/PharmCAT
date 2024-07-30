# PharmCAT

This is a customized version of PharmCAT specifically for use by SHC.

For standard PharmCAT, see https://github.com/PharmGKB/PharmCAT. 


## Building SHC PharmCAT

1. Run `make updateData`
2. Run `make updateShc`


## Update from standard PharmCAT

1. Update either `main` or `development` branch to match upstream PharmCAT branch
2. Rebase `shc_core` on top of this branch
3. Update SHC inputs as appropriate
4. Run `make updateData`
5. Run `make updateShc`
6. This becomes the new release


## Contact

For technical questions or bug reports, [file an issue](https://github.com/PharmGKB/PharmCAT/issues).

For general questions about the PharmCAT project, contact [pharmcat@pharmgkb.org](mailto:pharmcat@pharmgkb.org).


## Liability

:warning: PhamCAT assumes no responsibility for any injury to person or damage to persons or property arising out of, or related to any use of PharmCAT, or for any errors or omissions. The user recognizes they are using PharmCAT at their own risk.
