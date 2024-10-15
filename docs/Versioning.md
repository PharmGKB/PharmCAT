---
title: Versioning & Releases
permalink: versioning/
nav_order: 7
---
# Versioning & Releases

## Releases

As our team continues to work on PharmCAT, we will be releasing new versions and making them available on GitHub.
The latest release can always be found on the [GitHub releases page](https://github.com/PharmGKB/PharmCAT/releases).

Each release will include notes about changes that will affect users.  These notes might get a little technical.
If you prefer a simpler version or only care about major changes, take a look at the [Changelog](/changelog/).

We _highly recommend_ you use the **latest available** release of PharmCAT.


### Subscribing to release notifications

To get notifications for new releases, we recommend using [GitHub's notification system](https://docs.github.com/github/managing-subscriptions-and-notifications-on-github/configuring-notifications#:~:text=%22Configuring%20notifications.%22-,Configuring%20your%20watch%20settings%20for%20an%20individual%20repository,to%20being%20notified%20about%20issues.).

As of 10/14/2024, you can do this:

1. Go to the https://github.com/PharmGKB/PharmCAT
2. Select the Watch drop-down menu in the upper-right corner
3. Click Custom
4. Check Releases
5. Click Apply


## Versioning

Versions are identified by a version number (e.g. v1.6.0). The numbering follows the rules of
[Semantic Versioning](https://semver.org/). Each release also includes documentation describing
what has changed since the previous release.  This documentation follows the 
[Conventional Commit](https://www.conventionalcommits.org/en/v1.0.0/) guidelines. Since PharmCAT uses git as its
version control system, you can examine each change by looking through the
[git commit log](https://github.com/PharmGKB/PharmCAT/commits/development).


### When do we release a new version?

New versions of PharmCAT are released when there are one or more substantial changes that are deemed ready by the 
PharmCAT team. A substantial change can be one of the following:

- a new feature
- an update from a data source (e.g. PharmVar, CPIC)
- a fix for a problem reported by a user
- an improvement for an existing feature

**There is no time-based schedule for releasing new versions of PharmCAT**. 
