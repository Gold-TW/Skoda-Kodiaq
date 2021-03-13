#include "window.hpp"

MainWindow::MainWindow(QWidget *parent) : QWidget(parent) {
  main_layout = new QStackedLayout;

  homeWindow = new HomeWindow(this);
  main_layout->addWidget(homeWindow);

  settingsWindow = new SettingsWindow(this);
  main_layout->addWidget(settingsWindow);

  onboardingWindow = new OnboardingWindow(this);
  main_layout->addWidget(onboardingWindow);

  main_layout->setMargin(0);
  setLayout(main_layout);
  QObject::connect(homeWindow, SIGNAL(openSettings()), this, SLOT(openSettings()));
  QObject::connect(homeWindow, SIGNAL(closeSettings()), this, SLOT(closeSettings()));
  QObject::connect(homeWindow, SIGNAL(offroadTransition(bool)), this, SLOT(offroadTransition(bool)));
  QObject::connect(settingsWindow, SIGNAL(closeSettings()), this, SLOT(closeSettings()));

  // start at onboarding
  main_layout->setCurrentWidget(onboardingWindow);
  QObject::connect(onboardingWindow, SIGNAL(onboardingDone()), this, SLOT(closeSettings()));
  onboardingWindow->updateActiveScreen();

  // no outline to prevent the focus rectangle
  setStyleSheet(R"(
    * {
      font-family: Inter;
      outline: none;
    }
  )");
}

void MainWindow::offroadTransition(bool offroad){
  if(!offroad){
    closeSettings();
  }
}

void MainWindow::openSettings() {
  main_layout->setCurrentWidget(settingsWindow);
}

void MainWindow::closeSettings() {
  main_layout->setCurrentWidget(homeWindow);
}

bool MainWindow::eventFilter(QObject *obj, QEvent *event){
  if (event->type() == QEvent::MouseButtonPress) {
    homeWindow->glWindow->wake();
  }
  return false;
}
