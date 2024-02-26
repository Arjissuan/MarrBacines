from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains

class DbaaspSourceDownload:
    def __init__(self):
        self.options = Options()
        self.service = Service(executable_path='/snap/bin/geckodriver')
        self.url = 'https://dbaasp.org/search'
        self.driver = webdriver.Firefox(options=self.options, service=self.service)
        self.driver.maximize_window()
        self.driver.get(url=self.url)

    def getdata(self):
        pass


options = Options()
service = Service(executable_path='/snap/bin/geckodriver')
driver = webdriver.Firefox(options=options, service=service)

driver.maximize_window()
url = 'https://dbaasp.org/search'
driver.get(url)
assert 'Source' in driver.page_source
driver.implicitly_wait(6)
action = ActionChains(driver)
source_element = driver.find_element(by=By.XPATH, value='/html/body/main/div/div[3]/span/span/form/div[22]/span/span[1]/span')
search_element = driver.find_element(by=By.XPATH, value='/html/body/main/div/div[3]/span/div[1]')
source_element.location_once_scrolled_into_view
source_element.click()

for i in range(3):
    action.send_keys(Keys.ARROW_DOWN).perform()
action.send_keys(Keys.ENTER).perform()
print(source_element.text)

driver.implicitly_wait(3)
search_element.click()
export_data = driver.find_element(by=By.XPATH, value='/html/body/main/div/div[3]/div/div[1]/div[2]/a[1]')
export_data.click()