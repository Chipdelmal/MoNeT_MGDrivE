import time
import imghdr
import smtplib
from os import path
from email.message import EmailMessage

def mailAlert(
        password='1234',
        subject='Sim done!',
        sender='dummy@gmail.com',
        receiver='dummy@gmail.com',
        content='Sim finished!!!!',
        plotPath=None
    ):
    # Setup email text --------------------------------------------------------
    msg = EmailMessage()
    (msg['Subject'], msg['From'], msg['To']) = (subject, sender, receiver)
    msg.set_content(content)
    # Load and atttach plot ---------------------------------------------------
    if plotPath is not None:
        imgPath = path.join(plotPath)
        with open(imgPath, 'rb') as f:
            (image_data, image_type, image_name) = (
                f.read(), imghdr.what(f.name), f.name
            )
        msg.add_attachment(
            image_data, 
            maintype='image', subtype='image_type', filename=image_name
        )
    # Send email --------------------------------------------------------------
    with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
        smtp.login(sender, password)
        smtp.send_message(msg)