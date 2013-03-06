# VMWare Ubuntu image

We provide an Ubuntu 12.10 virtual machine image which runs under VMWare. To use this you will need at least 5.0GB of free space.

* Download and install [VMware Player](http://www.vmware.com/uk/products/desktop_virtualization/player/overview.html) 5.0.1 or above.
* Download our [virtual machine image](http://www2.epcc.ed.ac.uk/~michaelj/SoftwareCarpentry/UbuntuVM.zip)
 * Size: 1946346207
 * MD5 sum: 3faf8beb4780eea3bd75d452cff9a87d 
* Unzip UbuntuVM.zip
* Start VMware
* Click Open a Virtual Machine
* Browse into the unzipped directory
* Select Ubuntu.vmx
* Click Open
* Click Play virtual machine
* If a "This virtual machine might have been moved or copied." message appears then Click I moved it
* If a "Product VMware Player is running with a disk file..." message appears then Click OK
* If a "Cannot connect the virtual device floppy0" message appears then Click Yes
* Log-in page should now appear
* Enter password `ubuntu`
* Double-click Terminal icon at the bottom left of the screen, just below the Floppy Disk icon
* In the terminal window enter `ls`
* You should see some files and a `README.TXT` with information about installed software packages

The VM's root username is `vm` and password is `ubuntu`.
